from setuptools import setup

setup(
    use_scm_version={"write_to": "protlib_designer/_version.py"},
    setup_requires=['setuptools_scm'],
    entry_points={
        "console_scripts": [
            "protlib-designer=scripts.run_protlib_designer:run_protlib_designer",
            "protlib-plm-scorer=scripts.run_plm_scorer:run_plm_scorer",
        ]
    }
)
