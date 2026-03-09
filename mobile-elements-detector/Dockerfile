FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md ./

RUN pip install --no-cache-dir biopython rich python-dotenv matplotlib

COPY src/ ./src/

RUN pip install --no-cache-dir --no-deps .

ENV ENTREZ_EMAIL=""

CMD ["mobile-elements-detector"]