FROM python:2

COPY . /app/
RUN pip install --no-cache-dir -r /app/requirements.txt

WORKDIR /app

CMD ["python", "webserver.py"]
