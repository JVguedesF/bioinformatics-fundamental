FROM python:3.10-slim-bookworm

WORKDIR /app

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    pymol \
    wget \
    gcc \
    g++ \
    make \
    vienna-rna \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml .
COPY src ./src

RUN pip install --no-cache-dir ".[rna]"

CMD ["python", "src/main.py"]