from setuptools import setup

setup(
    name="ultimate-genomic-bot",
    version="1.0.0",
    py_modules=["bot", "init_db", "precision_medicine_pipeline"],
    entry_points={
        "console_scripts": [
            "genomic-bot=bot:main",
        ],
    },
    install_requires=[
        "cyvcf2>=0.30.0",
        "requests>=2.28.0",
    ],
    python_requires=">=3.8",
    author="Precision Medicine Engineering",
    description="A fully automated genomic and pharmacodynamic pipeline CLI tool.",
)