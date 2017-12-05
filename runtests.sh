if [[ -z $1 ]]; then
    flake8 bgwtools tests &&
    coverage run -m py.test tests &&
    coverage report
else
    flake8 bgwtools/$1.py tests/test_$1.py &&
    pytest -vv tests/test_$1.py
fi
