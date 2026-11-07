FROM python:3.13-slim-bookworm
RUN pip install --no-cache-dir numpy==2.1.3 scipy==1.14.1 matplotlib==3.9.2
ENV OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
COPY data/ /app/data/
RUN mkdir -p /app/output
WORKDIR /app
