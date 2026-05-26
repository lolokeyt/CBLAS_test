FROM ubuntu:22.04

# Установка необходимых зависимостей
RUN apt-get update && apt-get install -y \
    cmake \
    build-essential \
    pkg-config \
    libopenblas-dev \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

# Установка рабочей директории
WORKDIR /app

# Копирование всех файлов проекта
COPY . .

# Создание папки build и сборка проекта
RUN mkdir -p build && \
    cd build && \
    cmake .. && \
    make

# Команда по умолчанию (запуск тестов)
CMD ["./build/test_cblas"]
