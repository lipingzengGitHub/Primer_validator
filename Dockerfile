
FROM python:3.9-slim

# Install BLAST
RUN apt-get update && apt-get install -y ncbi-blast+ && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
RUN pip install biopython pandas openpyxl

WORKDIR /app
COPY primer_validator.py .
COPY example_data ./example_data

CMD ["python", "primer_validator.py", "--help"]



