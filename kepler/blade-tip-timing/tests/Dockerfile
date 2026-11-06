FROM python:3.13-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends curl && rm -rf /var/lib/apt/lists/*
RUN pip install pytest==8.4.1 pytest-json-ctrf==0.3.5 numpy==2.3.3 scipy==1.16.2 pandas==2.3.3
COPY . /tests/
# Ground truth readable by root only; the program under test runs as nobody.
RUN chown -R root:root /tests && chmod -R a+rX /tests/heldout && chmod 700 /tests/truth \
    && chmod 600 /tests/truth/* && mkdir -p /app
