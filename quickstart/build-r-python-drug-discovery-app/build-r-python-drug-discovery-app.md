author: Chanin Nantasenamat
id: build-r-python-drug-discovery-app
categories: snowflake-site:taxonomy/solution-center/certification/quickstart, snowflake-site:taxonomy/product/ai, snowflake-site:taxonomy/product/applications-and-collaboration, snowflake-site:taxonomy/snowflake-feature/snowpark, snowflake-site:taxonomy/industry/healthcare-and-life-sciences
language: en
summary: Train an R QSAR model for aromatase inhibition, serve it with R Plumber on Snowpark Container Services, and expose it through Streamlit and MCP.
environments: web
status: Published
feedback link: https://github.com/Snowflake-Labs/sfquickstarts/issues
fork repo link: https://github.com/sfc-gh-cnantasenamat/qsar-r-py-snowflake

# Build an R and Python Powered Drug Discovery App on Snowflake

## Overview

Many scientific models start life as R scripts on a laptop. This guide shows one path from that starting point to a shared, governed service on Snowflake without rewriting the science in another language.

You will work with AromatasePred, a QSAR (quantitative structure-activity relationship) model that predicts how strongly a molecule inhibits human aromatase (CYP19A1, ChEMBL target CHEMBL1978). Aromatase inhibitors such as letrozole and anastrozole are used to treat hormone-receptor-positive breast cancer.

R stays first-class for the science: `rcdk` computes molecular fingerprints, `ranger` trains a quantile random forest, and an R Plumber API serves the model in a container on Snowpark Container Services (SPCS). Python handles data access and the interface: Snowflake notebooks curate and explore ChEMBL data, a Streamlit in Snowflake app serves scientists, and an MCP (Model Context Protocol) endpoint serves AI assistants.

### What You'll Learn
- How to curate and explore ChEMBL bioactivity data in a Snowflake notebook
- How to featurize molecules with ECFP4 fingerprints and train a quantile random forest in R
- How to package an R model and Plumber API as a container image for SPCS
- How to call a public SPCS endpoint programmatically with a programmatic access token (PAT)
- How to front an R service with a Python Streamlit app and an MCP endpoint
- Which model limitations to disclose before sharing predictions

### What You'll Build
An end-to-end aromatase inhibitor predictor:

- A Snowflake notebook that turns raw ChEMBL records into 2,751 curated training compounds, with exploratory analysis
- An R training pipeline over those compounds
- An R Plumber service on SPCS with `/health`, `/predict`, and `/mcp` endpoints
- A multipage Streamlit app with Dashboard, Single Prediction, Batch Prediction, and About pages
- An MCP server exposing six tools to AI assistants

### What You'll Need
- Access to a [Snowflake account](https://signup.snowflake.com/?utm_source=snowflake-devrel&utm_medium=developer-guides&utm_cta=developer-guides) with ACCOUNTADMIN access for the one-time setup
- [Docker](https://www.docker.com/) with `linux/amd64` build support
- [Snowflake CLI](https://docs.snowflake.com/en/developer-guide/snowflake-cli/index) configured with a connection
- Python 3.11 or later with `requests`, `pandas`, `toml`, `rdkit`, and `snowflake-connector-python`
- Basic familiarity with R and SQL

## Setup

### Clone the Repository

```bash
git clone https://github.com/sfc-gh-cnantasenamat/qsar-r-py-snowflake.git
cd qsar-r-py-snowflake
```

The repository contains:

| Folder | Contents |
|---|---|
| `data/` | Snapshot of the cleaned ChEMBL1978 dataset used for the bundled model |
| `docker/` | R training script, Plumber API, Dockerfile, SPCS spec, bundled model artifacts |
| `streamlit/` | Multipage Streamlit app |
| `notebooks/` | `01_prepare_training_data.ipynb` (curation and EDA) and `02_score_compounds.ipynb` (dashboard predictions) |
| `native_app/` | Optional Native App package for per-consumer MCP deployment |
| `scripts/` | Endpoint smoke test and dashboard data generator |
| `sql/` | Snowflake setup, shared functions, optional agent tables and tools |
| `agent/` | Optional Cortex Agent specification and evaluation fixtures |

### Create Snowflake Objects

Open `sql/setup.sql` in a Snowsight SQL worksheet, replace `<your_user>` with your user name, and run steps 1 to 4. The script creates:

- `QSAR_ROLE`, a dedicated role that owns the project objects
- `QSAR_WH`, an extra-small warehouse
- `QSAR_DB` with a stage, an image repository (`QSAR_REPO`), and an audit table (`QSAR_PREDICTIONS_LOG`)
- `QSAR_CPU_POOL`, a single-node `CPU_X64_S` compute pool for the R service
- `QSAR_ST_POOL`, a small `CPU_X64_XS` pool for the Streamlit container runtime
- `QSAR_APP_EAI`, an external access integration that lets the app install packages from PyPI and call PubChem
- `QSAR_NB_EAI`, an external access integration that lets the notebooks call the ChEMBL API and install RDKit

```sql
USE ROLE ACCOUNTADMIN;

CREATE ROLE IF NOT EXISTS QSAR_ROLE;
GRANT ROLE QSAR_ROLE TO USER <your_user>;
-- SYSADMIN (and so ACCOUNTADMIN) inherits QSAR_ROLE, so admin CLI connections can use its objects
GRANT ROLE QSAR_ROLE TO ROLE SYSADMIN;

CREATE COMPUTE POOL IF NOT EXISTS QSAR_CPU_POOL
    MIN_NODES = 1
    MAX_NODES = 1
    INSTANCE_FAMILY = CPU_X64_S
    AUTO_RESUME = TRUE
    AUTO_SUSPEND_SECS = 3600
    COMMENT = 'QSAR aromatase SPCS inference';
GRANT USAGE, MONITOR ON COMPUTE POOL QSAR_CPU_POOL TO ROLE QSAR_ROLE;
GRANT BIND SERVICE ENDPOINT ON ACCOUNT TO ROLE QSAR_ROLE;
GRANT DATABASE ROLE SNOWFLAKE.PYPI_REPOSITORY_USER TO ROLE QSAR_ROLE;
```

`QSAR_ROLE` owns everything the guide creates. The SQL blocks in later steps start with `USE ROLE QSAR_ROLE`, so run them in a worksheet as written. Because `setup.sql` grants `QSAR_ROLE` to `SYSADMIN`, the Snowflake CLI commands also work from a connection that uses `SYSADMIN` or `ACCOUNTADMIN`, so you don't need a separate CLI connection.

**Check:** confirm that `QSAR_ROLE` can use the warehouse, both compute pools, both integrations, and the endpoint privilege:

```sql
USE ROLE QSAR_ROLE;
USE WAREHOUSE QSAR_WH;
SHOW GRANTS TO ROLE QSAR_ROLE;
SELECT "privilege", "granted_on", "name"
FROM TABLE(RESULT_SCAN(LAST_QUERY_ID()))
WHERE "granted_on" IN ('WAREHOUSE', 'COMPUTE_POOL', 'INTEGRATION', 'ACCOUNT')
ORDER BY 2, 3, 1;
```

You should see eight rows in this filtered output: `BIND SERVICE ENDPOINT` on the account, `USAGE` and `MONITOR` on `QSAR_CPU_POOL` and on `QSAR_ST_POOL`, `USAGE` on `QSAR_APP_EAI` and on `QSAR_NB_EAI`, and `USAGE` on `QSAR_WH`. The `SNOWFLAKE.PYPI_REPOSITORY_USER` database role is a separate grant that lets `RO5()` pin the notebook's RDKit version. If a row is missing, rerun the matching part of `setup.sql`.

Note the `repository_url` returned by `SHOW IMAGE REPOSITORIES`. You will push the container image there.

## Prepare Training Data

### Run the Notebook

`notebooks/01_prepare_training_data.ipynb` turns raw ChEMBL records into the training set, entirely inside Snowflake:

1. Fetches raw IC50 records for `CHEMBL1978` from the ChEMBL REST API into the `CHEMBL1978_RAW_ACTIVITIES` table: binding assays with an exact (`=`) relation
2. Keeps nM records with a SMILES string and a pChEMBL value, then keeps the highest pIC50 per molecule, saving the result as `CHEMBL1978_COMPOUNDS`
3. Labels a compound Active when pIC50 is at least 6.5, which corresponds to an IC50 of about 316 nM or lower
4. Explores the data: the pIC50 distribution, class balance, and Lipinski Rule of Five properties computed with RDKit
5. Writes `chembl1978_cleaned.csv` for training, plus `dashboard_eda.csv` and `dashboard_ro5.csv` for the app, to `@QSAR_STAGE`

Run it non-interactively from a SQL worksheet. Upload both notebooks to a stage, create a notebook project from that stage, then execute the first notebook on the app's compute pool with the `QSAR_NB_EAI` integration that `setup.sql` created:

```sql
USE ROLE QSAR_ROLE;
USE WAREHOUSE QSAR_WH;
USE SCHEMA QSAR_DB.PUBLIC;
CREATE STAGE IF NOT EXISTS NB_PROJECT_STAGE;
-- Upload the notebooks, for example:
-- snow sql -q "PUT file://notebooks/*.ipynb @QSAR_DB.PUBLIC.NB_PROJECT_STAGE AUTO_COMPRESS=FALSE OVERWRITE=TRUE"

CREATE OR REPLACE NOTEBOOK PROJECT QSAR_NOTEBOOKS FROM '@NB_PROJECT_STAGE';

EXECUTE NOTEBOOK PROJECT QSAR_NOTEBOOKS
    MAIN_FILE = '01_prepare_training_data.ipynb'
    COMPUTE_POOL = 'QSAR_ST_POOL'
    QUERY_WAREHOUSE = 'QSAR_WH'
    RUNTIME = 'V2-CPU-PY3.12'
    EXTERNAL_ACCESS_INTEGRATIONS = (QSAR_NB_EAI);
```

The run takes about three minutes. You can also open the notebook in a Snowsight Workspace and run it cell by cell on a Container Runtime with `QSAR_NB_EAI` attached, which displays the exploratory charts. If you used a database other than `QSAR_DB`, change `DATABASE` in the notebook's first code cell.

### Check the Result

From ChEMBL 37, the notebook keeps 2,751 compounds out of 3,960 raw records, 1,061 of them Active:

```sql
USE ROLE QSAR_ROLE;
USE WAREHOUSE QSAR_WH;
SELECT COUNT(*) AS N_COMPOUNDS, SUM(ACTIVITY_CLASS) AS N_ACTIVE
FROM QSAR_DB.PUBLIC.CHEMBL1978_COMPOUNDS;
```

A newer ChEMBL release can change these counts.

## Train the R Model

### Get the Training Data

Download the notebook's output into the repository's `data/` folder, replacing the snapshot there:

```bash
snow sql -q "GET @QSAR_DB.PUBLIC.QSAR_STAGE/chembl1978_cleaned.csv file://data/"
```

Each row is one compound:

```text
molecule_chembl_id,canonical_smiles,standard_value,pchembl_value,activity_class
CHEMBL4632445,CCC#CCOc1ccc2c(=O)c3ccc(Cn4ccnc4)cc3oc2c1,0.01,11.0,1
```

### Compute Fingerprints

`docker/train_rcdk.R` parses each SMILES with the Chemistry Development Kit (CDK) through `rcdk` and computes a 1024-bit ECFP4 circular fingerprint. Training and serving both use the same `rcdk` call, so features match exactly at inference time.

```r
fps <- lapply(mols, function(mol) {
  fp <- get.fingerprint(mol, type = "circular", circular.type = "ECFP4",
                        depth = 2, size = 1024)
  v  <- integer(1024)
  v[fp@bits] <- 1L
  v
})
```

Near-zero-variance bits are then removed with `caret::nearZeroVar()`, leaving 178 informative bits.

### Train a Quantile Forest

The script splits the data 80/20 and trains a 500-tree `ranger` random forest with `quantreg = TRUE`, which lets the service return prediction quantiles as well as a point estimate.

```r
model <- ranger(
  target ~ .,
  data       = train_df,
  num.trees  = 500,
  importance = "impurity",
  quantreg   = TRUE,
  seed       = 42
)
```

### Run Training in Docker

Build the serving image locally, then reuse it as the R runtime for training. This keeps the `rcdk`, Java, and `ranger` versions identical between training and serving.

```bash
docker build --platform linux/amd64 -f docker/Dockerfile.plumber -t r-qsar-plumber:local docker
mkdir -p output
docker run --rm --platform linux/amd64 \
  -v "$PWD/data:/data" -v "$PWD/output:/output" -v "$PWD/docker:/work" \
  r-qsar-plumber:local Rscript /work/train_rcdk.R
```

Training takes a few minutes under `linux/amd64` emulation. The script writes the model, the list of kept bits, the top 20 features, `train_fingerprints.csv`, and `training_metrics.json` to `output/`:

```json
{
  "oob_r2": 0.6588,
  "test_r2": 0.6537,
  "test_rmse": 0.7949,
  "n_compounds": 2751,
  "n_active": 1061,
  "n_train": 2203,
  "n_test": 548,
  "n_bits": 178,
  "n_trees": 500,
  "fingerprint": "CDK ECFP4 (rcdk 3.8.2)",
  "model": "ranger quantile RF"
}
```

Copy the artifacts into `docker/`, where the image build picks them up:

```bash
cp output/qsar_aromatase_model_rcdk.rds docker/qsar_aromatase_model.rds
cp output/kept_bits_rcdk.csv docker/kept_bits.csv
cp output/top_features_rcdk.csv output/training_metrics.json output/train_fingerprints.csv docker/
```

The repository also ships a bundled model, trained on the snapshot in `data/`, so you can skip this section and the previous one. That model reports test R2 0.6855 and RMSE 0.7708. The two snapshots contain the same compounds in a different row order, and the script draws its split by row position, so the metrics differ. Read the Review Model Limitations section before quoting either set of numbers.

## Serve With R Plumber

### Inspect the API

`docker/plumber_api.R` loads the model and measured-compound fingerprints once at startup. Its core endpoints are:

| Endpoint | Purpose |
|---|---|
| `GET /health` | Readiness probe; returns tree and feature counts |
| `POST /predict` | Batch prediction in the Snowflake service-function row format |
| `POST /mcp` | MCP server over JSON-RPC 2.0 |

Each prediction returns the median and the 5th and 95th percentiles from the quantile forest:

```r
pred   <- predict(model, data = fp_df, type = "quantiles",
                  quantiles = c(0.05, 0.5, 0.95))
```

### Build and Push

The image uses `rocker/r-ver:4.3.2` with a JDK for `rcdk` and bakes in the 4 MB model, so no stage mount is needed at runtime.

```bash
cd docker
snow spcs image-registry login
docker build --platform linux/amd64 -f Dockerfile.plumber \
  -t <repository_url>/r-qsar-plumber:v1 .
docker push <repository_url>/r-qsar-plumber:v1
cd ..
```

### Create the Service

Create the service from `docker/plumber_spec.yaml`:

```sql
USE ROLE QSAR_ROLE;
USE SCHEMA QSAR_DB.PUBLIC;

CREATE SERVICE QSAR_PLUMBER_SVC
    IN COMPUTE POOL QSAR_CPU_POOL
    FROM SPECIFICATION $$
spec:
  containers:
    - name: qsar-r-plumber
      image: /qsar_db/public/qsar_repo/r-qsar-plumber:v1
      env:
        PORT: "5000"
      readinessProbe:
        port: 5000
        path: /health
      resources:
        requests:
          memory: 0.5Gi
          cpu: 0.5
        limits:
          memory: 13Gi
          cpu: 3
  endpoints:
    - name: api
      port: 5000
      public: true
    $$
    MIN_INSTANCES = 1
    MAX_INSTANCES = 1;
```

Wait for the container to report `READY`, then get the public URL:

```sql
USE ROLE QSAR_ROLE;
SHOW SERVICE CONTAINERS IN SERVICE QSAR_DB.PUBLIC.QSAR_PLUMBER_SVC;
SHOW ENDPOINTS IN SERVICE QSAR_DB.PUBLIC.QSAR_PLUMBER_SVC;
```

The container logs should show the model loading:

```sql
USE ROLE QSAR_ROLE;
USE WAREHOUSE QSAR_WH;
SELECT SYSTEM$GET_SERVICE_LOGS('QSAR_DB.PUBLIC.QSAR_PLUMBER_SVC', 0, 'qsar-r-plumber', 20);
```

```text
Loading model artifacts...
  ranger model: 500 trees, 178 features
  training fingerprints: 2751 compounds
Warming up rcdk JVM...
  rcdk JVM ready
Plumber API ready — listening on port 5000
```

### Call the Endpoint

Public SPCS endpoints require Snowflake authentication. For programmatic clients, generate a PAT in Snowsight under your user profile, then pass it in the `Authorization` header in this exact form:

```text
Authorization: Snowflake Token="<your-PAT>"
```

Sending the PAT as `Bearer <your-PAT>` returns a redirect to the browser sign-in page instead.

`/predict` expects rows of `[row_index, smiles]` under a `data` key. A body such as `{"smiles": "..."}` returns `{"data": []}` with HTTP 200, so check the request shape if you get empty results.

Run the included smoke test:

```bash
export QSAR_ENDPOINT="<ingress_url from SHOW ENDPOINTS>"
export SNOWFLAKE_PAT="<your-PAT>"
python scripts/test_endpoint.py
```

Example output:

```text
health: 200 {"status":"ok","trees":500,"features":178}

compound      pIC50  class        q05   q95
letrozole      8.26  Active      5.79  9.22
anastrozole    8.64  Active      6.95  9.53
exemestane     6.21  Inactive    4.96  7.83
ethanol        5.21  Inactive    4.22  6.92
invalid           -  Invalid SMILES: not_a_smiles     -     -
```

## Build The Streamlit Frontend

### Review the App

The `streamlit/` folder contains a multipage app. It runs in Streamlit in Snowflake, on Streamlit Community Cloud, or locally:

- **Dashboard** charts the training data, Lipinski Rule of Five compliance, and model predictions
- **Single Prediction** predicts pIC50 for one molecule, entered as a SMILES string or picked from example molecules, and logs the result to `QSAR_PREDICTIONS_LOG`
- **Batch Prediction** scores many SMILES at once, starting from the same example molecules
- **About** documents the architecture

`streamlit/_shared.py` calls the R service over the SPCS internal DNS name, using the same batch payload as the smoke test:

```python
resp = requests.post(
    get_predict_url(),
    headers=headers,
    json={"data": [[0, smiles]]},
    timeout=60,
)
```

### Check the App Service

The app reads the service DNS name using `SHOW SERVICES` when it runs in Snowflake, so you do not need to edit a URL. With the default names, `streamlit/_shared.py` is ready to deploy. If you use a different database or schema, update `_DB_SCHEMA` in that file. Local development requires `spcs.predict_url` in Streamlit secrets.


### Score Compounds

The Dashboard compares measured and predicted pIC50 for every compound. First, run `sql/functions.sql` in Snowsight after the service is ready. The file creates `PREDICT_PIC50`, `ECFP4`, `SIMILAR_COMPOUNDS`, and `EXPLAIN_PREDICTION` over the R service, plus an RDKit `RO5` function. The setup grants `QSAR_ROLE` access to Snowflake's PyPI repository for the same RDKit version used in the first notebook. Select **Check** at the end of `sql/functions.sql`: the letrozole prediction should be about 7.89 with the bundled model. If you trained on notebook output, your result may differ.

`notebooks/02_score_compounds.ipynb` then reads `CHEMBL1978_COMPOUNDS` and scores every compound through the SQL functions, which call the R service internally without a token. It creates predictions, fingerprints, and descriptors tables and writes `dashboard_predictions.csv` to `@QSAR_STAGE`. Run it from the same notebook project:

```sql
USE ROLE QSAR_ROLE;
USE WAREHOUSE QSAR_WH;
EXECUTE NOTEBOOK PROJECT QSAR_DB.PUBLIC.QSAR_NOTEBOOKS
    MAIN_FILE = '02_score_compounds.ipynb'
    COMPUTE_POOL = 'QSAR_ST_POOL'
    QUERY_WAREHOUSE = 'QSAR_WH'
    RUNTIME = 'V2-CPU-PY3.12'
    EXTERNAL_ACCESS_INTEGRATIONS = (QSAR_NB_EAI);
```

Then upload the feature importances and metrics that the About and Dashboard pages read:

```bash
cp docker/top_features_rcdk.csv /tmp/top_features.csv
snow sql -q "PUT file:///tmp/top_features.csv @QSAR_DB.PUBLIC.QSAR_STAGE AUTO_COMPRESS=FALSE OVERWRITE=TRUE"
snow sql -q "PUT file://docker/training_metrics.json @QSAR_DB.PUBLIC.QSAR_STAGE AUTO_COMPRESS=FALSE OVERWRITE=TRUE"
```

If you skipped the notebooks and used the bundled model, generate all three dashboard files from your laptop instead. `scripts/generate_dashboard_data.py` reads `data/chembl1978_cleaned.csv` and calls the service's public endpoint with a PAT:

```bash
export QSAR_ENDPOINT="<ingress_url>"
export SNOWFLAKE_PAT="<your-PAT>"
export SNOWFLAKE_CONNECTION="<your connection name>"
python scripts/generate_dashboard_data.py
```

Set `QSAR_STAGE` as well if your stage is not `@QSAR_DB.PUBLIC.QSAR_STAGE`.

### Deploy

The app installs `rdkit` and `streamlit-jsme` from PyPI, so it runs on the Streamlit container runtime. `streamlit/snowflake.yml` already declares the runtime, compute pool, and external access integration:

```yaml
    runtime_name: SYSTEM$ST_CONTAINER_RUNTIME_PY3_11
    compute_pool: QSAR_ST_POOL
    external_access_integrations:
      - QSAR_APP_EAI
```

From the `streamlit/` folder, deploy with the Snowflake CLI:

```bash
cd streamlit
snow streamlit deploy --replace
```

Open the app from the URL the command prints. The first load takes a minute or two while the container installs packages. On Single Prediction, the **Input** panel offers example molecules: three FDA-approved aromatase inhibitors (letrozole, anastrozole, and exemestane) and two inactive compounds (aspirin and caffeine). Pick one, or paste your own SMILES string, then select **Predict pIC50**. The result shows the predicted pIC50, its quantile range, and the activity class, and each prediction is written to `QSAR_PREDICTIONS_LOG`.

The page also includes a JSME structure editor for drawing molecules. It appears only when the app runs locally or on Community Cloud. In Snowflake, `running_in_snowflake()` in `_shared.py` detects the active Snowpark session and the page shows SMILES input only.

After a redeploy, the running app can keep serving the previous code. If changes do not appear, restart the app from Snowsight.

The app's PubChem name and fallback-image lookups are off by default. If you opt in on Single Prediction, a SMILES-derived InChIKey and possibly the SMILES itself leave Snowflake; keep this unchecked for proprietary structures. In Streamlit in Snowflake these calls also need the external access integration.

## Connect AI Assistants

The same container exposes an MCP server at `/mcp`, so AI assistants can call the R model directly. It offers six tools:

| Tool | Purpose |
|---|---|
| `predict_pic50` | Predict pIC50, activity class, and quantile range for a SMILES string |
| `validate_smiles` | Check a SMILES string with the CDK parser |
| `get_model_info` | Return model metadata and recorded metrics |
| `get_dataset_stats` | Summarize the training dataset |
| `get_feature_importance` | Rank fingerprint bits by impurity importance |
| `get_ro5_stats` | Report Rule of Five compliance by activity class |

Register the server in an MCP client that supports HTTP transport. For Cortex Code, add an entry to `~/.snowflake/cortex/mcp.json`:

```json
{
  "mcpServers": {
    "qsar-aromatase": {
      "transport": "http",
      "url": "https://<your-ingress-url>/mcp",
      "headers": {
        "Authorization": "Snowflake Token=\"<your-PAT>\""
      }
    }
  }
}
```

Then ask questions such as "Predict the aromatase activity of letrozole, SMILES `N#Cc1ccc(cc1)C(n1cncn1)c1ccc(C#N)cc1`". A `predict_pic50` call returns:

```text
SMILES: N#Cc1ccc(cc1)C(n1cncn1)c1ccc(C#N)cc1
pIC50:  8.26  (Active)
90% CI: [5.79, 9.22]
IC50 ~  5.5 nM
```

`docs/mcp-walkthrough.md` in the repository covers more clients and an optional Native App package that lets each consumer run a private copy of the service.

## Ask an Agent

This optional section describes a prototype tested in isolated `QSAR_AG_*` objects. It is not deployed to the live `QSAR_DB` app. After completing the core flow, you can add an in-platform Cortex Agent to your own guide environment. Notebook 01 also writes `CHEMBL1978_ASSAYS` and `CHEMBL1978_DOCUMENTS` from ChEMBL. In Snowsight, run `sql/agent.sql` to create a semantic view, a Cortex Search service over linked assay descriptions and publication abstracts, and a procedure that ranks up to 20 candidate SMILES. The search service adds indexing costs.

Create the agent from `agent/aromatase_agent.yaml` using the companion repository's deployment script, with the Snowflake connection you used for the core flow:

```bash
python scripts/deploy_agent.py --connection <your_connection>
```

The script creates `QSAR_DB.PUBLIC.AROMATASE_AGENT` and adds it to Snowflake Intelligence. Restart or redeploy the app to see **Ask an Agent**. For a first question, ask: "How many measured aromatase compounds are active and inactive?" The expected answer is 1,061 active and 1,690 inactive. Then ask about the assay conditions for exemestane. Search should cite ChEMBL assay identifiers.

The agent uses Analyst for measured data, Search for assay and publication text, and R-backed tools for predictions and explanations. Its similarity calculation compares against all measured compounds, including held-out rows, and does not establish a validated applicability domain. The fingerprint explanation ranks globally important bits present in a molecule; it is not a per-molecule attribution. `agent/eval_questions.sql` and `agent/eval_config.yaml` provide an optional evaluation set. Its expected prediction values use the bundled model; update them if you trained on notebook output. In the isolated prototype, all 15 cases completed with mean answer correctness 0.956, logical consistency 0.956, and tool selection 0.967. Those scores describe agent responses and routing, not model potency accuracy.

## Review Model Limitations

This is a teaching example, not a validated drug discovery model. Disclose these limitations wherever you share predictions:

- **Feature filtering before the split.** `nearZeroVar()` runs on all compounds before the train/test split, so the test set influences which bits are kept. Move the filter after the split and fit it on training data only for a leakage-safe estimate.
- **One split, order-dependent.** The 80/20 split is drawn by row position with a fixed seed. The same compounds in a different row order give test R2 0.6537 instead of 0.6855, so treat a single split's metrics as approximate and prefer cross-validation.
- **Reported R2 is squared correlation.** `test_r2` is computed as `cor(pred, y)^2`, not the conventional `1 - SSE/SST`. Report both, or use the conventional definition.
- **Quantiles are not calibrated intervals.** The API labels the 5th to 95th percentile range a "90% CI", but it has not been checked for coverage on held-out data.
- **No validated applicability domain.** Ethanol receives a model prediction rather than a validated out-of-domain warning. The optional agent reports nearest-measured-compound similarity as a heuristic, not an applicability-domain test.
- **Held-out misses on known drugs.** Exemestane, a steroidal inhibitor with a measured pIC50 of 7.82 in the dataset, falls in the notebook-trained model's held-out test set and is predicted at 6.21, so it is labelled Inactive. The bundled model's training split includes exemestane and predicts 7.00. Steroids are common in the data (682 of 2,751 compounds) but mostly inactive (31% active, against 39% overall), so check predictions for potent steroids against measured values.
- **Audit identity.** In Streamlit in Snowflake, `USER_NAME` in `QSAR_PREDICTIONS_LOG` records the app's service user, not the person making the prediction. Capture the viewer's identity explicitly before relying on the table for auditing.
- **Silent empty responses.** A malformed `/predict` request returns HTTP 200 with no rows. Validate the request schema before relying on the service in automated pipelines.

## Clean Up

Compute pools bill while they have active nodes. Remove the objects created in this guide:

If you completed the optional Agent section, first remove it from Snowflake Intelligence and drop the agent and search service:

```sql
USE ROLE QSAR_ROLE;
ALTER SNOWFLAKE INTELLIGENCE SNOWFLAKE_INTELLIGENCE_OBJECT_DEFAULT
  DROP AGENT QSAR_DB.PUBLIC.AROMATASE_AGENT;
DROP AGENT IF EXISTS QSAR_DB.PUBLIC.AROMATASE_AGENT;
DROP CORTEX SEARCH SERVICE IF EXISTS QSAR_DB.PUBLIC.ASSAY_SEARCH;
```

Then clean up the core project:

```sql
USE ROLE QSAR_ROLE;
DROP STREAMLIT IF EXISTS QSAR_DB.PUBLIC.QSAR_AROMATASE_APP;
DROP SERVICE IF EXISTS QSAR_DB.PUBLIC.QSAR_PLUMBER_SVC;

USE ROLE ACCOUNTADMIN;
DROP COMPUTE POOL IF EXISTS QSAR_CPU_POOL;
DROP COMPUTE POOL IF EXISTS QSAR_ST_POOL;
DROP INTEGRATION IF EXISTS QSAR_APP_EAI;
DROP INTEGRATION IF EXISTS QSAR_NB_EAI;
DROP DATABASE IF EXISTS QSAR_DB;
DROP WAREHOUSE IF EXISTS QSAR_WH;
DROP ROLE IF EXISTS QSAR_ROLE;
```

Also revoke any PAT you created for this guide from your user profile in Snowsight.

## Conclusion And Resources

Congratulations! You've successfully trained an R QSAR model, served it with R Plumber on Snowpark Container Services, and exposed it to scientists through Streamlit and to AI assistants through MCP.

### What You Learned
- How to keep R as the scientific engine while Python provides the interface
- How to package an R model and API as an SPCS service
- How to authenticate to a public SPCS endpoint with a PAT
- How one container can serve both an app and an MCP server
- Which validation gaps to close before treating a model as production-ready

### Related Resources
- [Companion repository](https://github.com/sfc-gh-cnantasenamat/qsar-r-py-snowflake)
- [Snowpark Container Services](https://docs.snowflake.com/en/developer-guide/snowpark-container-services/overview)
- [Access a public endpoint programmatically](https://docs.snowflake.com/en/developer-guide/snowpark-container-services/tutorials/advanced/tutorial-8-access-public-endpoint-programmatically)
- [Programmatic access tokens](https://docs.snowflake.com/en/user-guide/programmatic-access-tokens)
- [Streamlit in Snowflake](https://docs.snowflake.com/en/developer-guide/streamlit/about-streamlit)
- [ChEMBL database](https://www.ebi.ac.uk/chembl/)
