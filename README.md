# pprint
## Setup

ppprint requires Python3 and uses the webframework [Django](https://www.djangoproject.com/start/).

Make sure you have the following dependencies installed:
- Redis

To install them on Ubuntu/Debian use:
```shell
apt-get install redis-server
```

### Other dependencies
Python dependencies are managed via pip. It is however advised to install those dependencies into an virtual environment. This can be created with:
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

## Usage

To set up the database, run `python manage.py migrate`.

Jobs are managed with [Celery](https://docs.celeryproject.org/en/stable/), a distributed Task-Queue for Python. Start a worker with:
```shell
celery -A halfpipe worker --task-events
```

The Django Webserver can be started with:
```shell
python manage.py runserver
```

### Best practices

The Python code style conforms to [Black](https://github.com/psf/black). To format all Python code use:
```shell
black .
```
