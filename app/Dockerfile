FROM python:3.13

WORKDIR /app
COPY . .

RUN apt-get update && apt-get install -y wget unzip

RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.16/IQ-TREE-2.3.16-Linux.tar.gz \
    && tar -xvzf IQ-TREE-2.3.16-Linux.tar.gz \
    && mv iqtree-2.3.16-Linux/bin/iqtree3 /usr/local/bin/ \
    && rm -rf iqtree-2.3.16-Linux IQ-TREE-2.3.16-Linux.tar.gz

RUN pip install --upgrade pip setuptools wheel
RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8080
ENV PORT=8080

CMD ["gunicorn", "-b", "0.0.0.0:8080", "main:app"]
