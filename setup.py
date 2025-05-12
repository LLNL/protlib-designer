from setuptools import setup

setup(
    use_scm_version={
        "write_to": "protlib_designer/_version.py",
        "local_scheme": "no-local-version",
    },
    setup_requires=['setuptools_scm']
)
