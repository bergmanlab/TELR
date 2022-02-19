import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TELR",
    version="1.0",
    author="Shunhua Han",
    author_email="hanshunhua0829@gmail.com",
    description="A fast non-reference transposable element detector from long read sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bergmanlab/TELR",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    entry_points={"console_scripts": ["telr = telr.telr:main"]},
    classifiers=(
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ),
)