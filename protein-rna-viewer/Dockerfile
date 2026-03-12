FROM python:3.10-slim-bookworm

WORKDIR /app

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget gcc g++ make libgsl-dev \
    pymol libgl1-mesa-glx libqt5gui5 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV QT_QPA_PLATFORM=offscreen
ENV DISPLAY=:0
ENV PYTHONPATH=/app

COPY pyproject.toml .
RUN pip install --no-cache-dir ".[rna]"

COPY src ./src

RUN mkdir -p /app/data /app/results

CMD ["python", "src/main.py"]