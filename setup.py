import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DNA-Toolkit",
    version="1.0.0",
    author="bobroxq",
    author_email="bobroxq@gmail.com",
    description="A small suite for DNA manipulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bobroxq/DNA_Toolkit",
    project_urls={
        "Bug Tracker": "https://github.com/bobroxq/DNA_Toolkit/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
)