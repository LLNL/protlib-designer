from setuptools import setup

setup(
    entry_points={
        "console_scripts": [
            "protlib-designer=scripts.run_protlib_designer:run_protlib_designer",
            "protlib-plm-scorer=scripts.run_plm_scorer:run_plm_scorer",
        ]
    }
)
