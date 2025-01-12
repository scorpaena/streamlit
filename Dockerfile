FROM python:3.10-slim-bullseye

WORKDIR /app
COPY . /app
RUN apt-get update && apt-get -y install < packages.txt

RUN pip install --no-cache-dir -r requirements.txt

ENV OS_TYPE=linux

EXPOSE 8501

ENTRYPOINT ["streamlit", "run", "gear_model.py", "--server.port=8501", "--server.address=0.0.0.0"]