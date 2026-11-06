# Workstation image for the rising-main surge study.
FROM python:3.13-slim-bookworm

ARG PIP_NO_CACHE_DIR=1
RUN python3 -m pip install numpy==2.2.6 scipy==1.15.3 pandas==2.2.3 matplotlib==3.10.3

ENV OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 MPLBACKEND=Agg
COPY data /app/data
WORKDIR /app
RUN mkdir -p /app/output /app/work
