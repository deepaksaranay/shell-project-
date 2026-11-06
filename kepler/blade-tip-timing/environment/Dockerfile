FROM python:3.13-slim

RUN pip install --no-cache-dir numpy==2.3.3 scipy==1.16.2 pandas==2.3.3 pytest==8.4.1

WORKDIR /app
COPY data/ /app/data/
RUN mkdir -p /app/output
