# load modules
import sys
import pickle
import numpy as np
import warnings
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler

# Get drug name from command-line argument
if len(sys.argv) != 3:
    print("Usage: python train_elastic_net.py <drug_name> <training_sets_dict.pkl>")
    sys.exit(1)

drug_name = sys.argv[1]
training_sets_file = sys.argv[2]
save_path = "/group/iorio/Raffaele/SCDRESP_data/data/results/elnet_nopoly_rep0/"
run_random=False

# Load the training sets dictionary
with open(training_sets_file, "rb") as f:
    training_sets_dict = pickle.load(f)

if drug_name not in training_sets_dict:
    print(f"Error: {drug_name} not found in training sets.")
    sys.exit(1)

data = training_sets_dict[drug_name]

print(f"Training ElasticNet for drug: {drug_name}...", flush=True)

try:
    warnings.simplefilter(action="ignore", category=FutureWarning)

    X, Y = data["X"], data["Y"]
    X = X.fillna(0)
    Y = Y.fillna(Y.mean()).to_numpy() 
    X_features = X.columns.to_numpy()

    if run_random:
        np.random.seed(42)  
        Y = np.random.permutation(Y)
     
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Remove highly correlated features
    cmat = np.corrcoef(X_scaled, rowvar=False)
    np.fill_diagonal(cmat, 0)
    idx_hcf = list(np.where(np.abs(cmat) > 0.99)[1])# changed because had issues of convergence
    X_scaled = np.delete(X_scaled, idx_hcf, axis=1)
    retained_X_features = np.delete(X_features, idx_hcf)
    
    # Define nested cross-validation
    outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
    inner_cv = KFold(n_splits=3, shuffle=True, random_state=42)

    # Hyperparameter tuning
    elnet = ElasticNetCV(
        alphas=np.logspace(-4, 1, 50),
        l1_ratio=np.linspace(0.1, 1.0, 10),
        max_iter=10000,
        tol=1e-4,
        cv=inner_cv,
        random_state=42,
        n_jobs=-1
    )

    # Outer cross-validation
    r2_scores = cross_val_score(elnet, X_scaled, Y, scoring='r2', cv=outer_cv, n_jobs=-1)
    rmse_scores = np.sqrt(-cross_val_score(elnet, X_scaled, Y, scoring='neg_mean_squared_error', cv=outer_cv, n_jobs=-1))

    # Train on full dataset
    elnet.fit(X_scaled, Y)
    best_alpha = elnet.alpha_
    best_l1_ratio = elnet.l1_ratio_

    # Save the model
    model_dict = {
        "model": elnet,
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "best_alpha": best_alpha,
        "best_l1_ratio": best_l1_ratio,
        "cv_r2_scores": r2_scores,
        "mean_r2": np.mean(r2_scores),
        "std_r2": np.std(r2_scores),
        "cv_rmse_scores": rmse_scores,
        "mean_rmse": np.mean(rmse_scores),
        "std_rmse": np.std(rmse_scores),
        "pre_filter_features":X_features,
        "post_filter_features":retained_X_features
    }

    model_filename = f"{save_path}/elastic_net_model_{drug_name}.pkl"
    with open(model_filename, "wb") as f:
        pickle.dump(model_dict, f)

    print(f"Saved model for {drug_name}: {model_filename}")

except Exception as e:
    print(f"Error training {drug_name}: {e}", flush=True)
    sys.exit(1)

