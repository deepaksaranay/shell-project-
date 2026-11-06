FROM python:3.13-slim-bookworm

ENV OMP_NUM_THREADS=2 \
    MKL_NUM_THREADS=2 \
    OPENBLAS_NUM_THREADS=2

RUN pip install --no-cache-dir numpy==2.3.3 scipy==1.16.2 pandas==2.3.3

WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output
