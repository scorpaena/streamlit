FROM python:3.10

COPY . /app
WORKDIR /app

RUN apt-get update && apt-get install -y libgl1 libxrender1 xvfb python3-vtk9 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

ENTRYPOINT ["streamlit", "run", "gear_model.py", "--server.port=8501", "--server.address=0.0.0.0"]