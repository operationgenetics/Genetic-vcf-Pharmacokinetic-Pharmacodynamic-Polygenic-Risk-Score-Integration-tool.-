from setuptools import setup

setup(
    name="ultimate-genomic-bot",
    version="1.0.0",
    py_modules=["bot"],
    entry_points={
        "console_scripts": [
            "genomic-bot=bot:main",
        ],
    },
    install_requires=[],
    author="Your Name",
    description="A fully automated genomic and pharmacodynamic pipeline CLI tool.",
)
