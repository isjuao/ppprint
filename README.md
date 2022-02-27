# pprint
## Setup
For development, create virtual environment with `python3 -m venv venv`. Later use `source venv/bin/activate`.
Install requirements with `pip install -r requirements/dev.txt`.

## Usage
### Local testing: website
For now, run Django project `ppprint` from PyCharm. Start celery worker in separate terminal with `celery -A ppprint worker`.

### Local testing: data only
Use `from ppprint.preprocessing import run` followed by `run(6)` to quickly test data handling of a subproteome of E. coli 
located in `media/import_job/6/`, without using the website.
