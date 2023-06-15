FROM python:3.9

WORKDIR /app

COPY requirements.txt /app

RUN pip install --no-cache-dir -r requirements.txt

COPY . /app

ENV PYTHONPATH "${PYTHONPATH}:/app"

CMD python refgenDetector/refgenDetector_main.py

https://stackoverflow.com/questions/46245844/pass-arguments-to-python-argparse-within-docker-container