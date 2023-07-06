# ppprint
**P**redict**P**rotein **print**

A web application that analyses [PredictProtein](https://predictprotein.org/) prediction data for whole proteomes of multiple organisms and 
displays 
the visualised comparison results in feature-specific dashboards. Protein features currently available include 
protein disorder, transmembrane helices, secondary structure and binding sites.

For background information, see below setup instructions.

## Setup and Usage

**ppprint** requires Python3 and uses the webframework [Django](https://www.djangoproject.com/start/).

Make sure you have Redis installed. To install it on Ubuntu/Debian use:
```shell
apt-get install redis-server
```

Python **dependencies** are managed via pip. It is however advised to install those dependencies into an virtual 
environment. This can be created with:
```shell
python3 -m venv venv
```

To activate the environment use:
```shell
source venv/bin/activate
```
This needs to be done in every shell that should use our newly created environment.

For development purposes we suggest installing the Python dependencies in the `requirements/dev.txt` file:
```shell
pip install -r requirements/dev.txt
```

To set up the database, run:

```shell
python manage.py migrate
```

Jobs are managed with [Celery](https://docs.celeryproject.org/en/stable/), a distributed Task-Queue for Python. Start a worker with:
```shell
celery -A halfpipe worker --task-events
```

The Django webserver can be started with:
```shell
python manage.py runserver
```

## Background

PredictProtein is a collection of a multitude of protein feature prediction tools. While
for collections of proteins such as whole proteomes, generating predictions requires a local
installation of PredictProtein, the service can be accessed online for single sequences.
Here, interactive dashboards present the results for several protein features including
those mentioned above. Few forms of visualizations display basic
information such as the location of the feature regions, for instance transmembrane helices,
within the concerned protein. Since predictions were produced only for the single sequence,
no distributions as for collections of proteins can be calculated. However, for multiple
sequences up to whole proteomes and thus using PredictProtein locally instead of the web
service, no additional statistics or feature distributions are provided. Here, **ppprint** serves as a potential 
extension for PredictProtein.

The aim of **ppprint** is to calculate trends and distributions among whole-proteome data. Analysis results as well 
as comparisons of up to four proteomes are then accessible via per-feature dashboards provided by the
web application. Users can upload their own data or use available sample proteomes to
examine underlying patterns by exploring the provided forms of visualization for protein
disorder, transmembrane helices and protein binding sites. Derived insights into how
single organisms or groups behave regarding crucial protein features can then enhance
knowledge on their relationship and thus aid in research and drug development. Because
of its non-final state, the **ppprint** website can currently only be loaded locally following
the instructions as given in the git repository. In later stages of development, the web
application will be accessible as an online service.

### Using a local version of ppprint

#### Uploading data

In order to make PredictProtein predictions available for comparison, the user has to
upload them via the upload page, which can be reached using the navigation
bar. Since for most proteomes of model organisms, Gigabytes of prediction files may have
been generated, the user can only upload a compressed .tar file with an extension such
as .xz, .gz. For a successful data import, the user also has to enter a name and may pick
a color for the display of analysis results. When the user has selected a file and filled out
the remaining two fields, the data can be submitted for internal parsing and extraction.

#### Comparing proteomes

The success of an upload can be verified via the proteome selection page.
When the specified name of the uploaded proteome appears in the list of proteomes
selectable for comparison, data import and extraction has finished successfully. The user
can then choose one or multiple proteomes and submit them for analysis. Currently, the
maximum number of proteomes that can be simultaneously selected is four, ensuring that
there is no visual overload in the final plot figures.

Once the user has submitted a set of proteomes for analysis, an overview of the created
comparisons can be accessed by navigating to the listing page. Here, the
creation date of the job and its status, which can be "created", "running", "success", and
"failure" are given for each created job. The user can then request the final visualization of
analysis results by clicking on the respective comparison identifier, internally corresponding
to the primary key of the VisualizationJob. If the status of the job is not "success", the
list of plots ready to be displayed may be incomplete and either requires a page refresh ("created"/"running") or starting a new comparison ("failure"), for analysis results to be
shown completely.

Finally, the results of the analysis performed by **ppprint** can be accessed via the detail
result pages. The overview tab displays feature-independent properties such
as protein lengths and sample sizes of the selected proteomes. Using the tab bar, the user can navigate to the 
feature-related detail pages in order to view the produced forms of visualization assembled as
basic dashboards. Since **ppprint** is still in development, user-guiding indications such as
plot-specific help texts for interpretation and the introduction of required concepts are yet
to be added.
