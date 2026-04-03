FROM python:3.13

RUN apt-get update && apt-get install -y --no-install-recommends \
    mafft \
    muscle \
    wget \
    tar \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y wget tar \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com \
    chmod +x muscle-linux-x86.v5.3 \
    ./muscle-linux-x86.v5.3

RUN wget https://github.com/iqtree/iqtree3/releases/download/v3.1.1/iqtree-3.1.1-Linux.tar.gz \
    && tar -xvzf iqtree-3.1.1-Linux.tar.gz \
    && mv iqtree-3.1.1-Linux/bin/* /usr/local/bin/ \
    && chmod +x /usr/local/bin/iqtree3

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PORT=8080

CMD gunicorn main:app -b 0.0.0.0:$PORT -w 2 --timeout 300
