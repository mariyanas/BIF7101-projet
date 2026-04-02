FROM python:3.11-slim

RUN apt-get update && apt-get install -y wget tar

RUN wget http://www.iqtree.org/doc/iqtree-2.2.2-Linux.tar.gz \
    && tar -xvzf iqtree-2.2.2-Linux.tar.gz \
    && mv iqtree-2.2.2-Linux/bin/iqtree3 /usr/local/bin/ \
    && chmod +x /usr/local/bin/iqtree3

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . /app
WORKDIR /app

CMD ["gunicorn", "main:app"]
