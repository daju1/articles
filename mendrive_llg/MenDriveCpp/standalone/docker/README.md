# Docker для MenDrive standalone (Ubuntu)

Образ на Ubuntu 22.04: g++, make, OpenMP (`libomp-dev`), `/usr/bin/time`, Python 3 с numpy, matplotlib, pandas, scipy, Jupyter.

## Установка Docker на Ubuntu (один раз)
```sh
sudo apt-get install -y docker.io
sudo usermod -aG docker $USER   # затем перелогиниться
```

## Запуск
```sh
cd MenDriveCpp/standalone/docker
chmod +x run.sh
./run.sh build                       # собрать образ (~2–3 мин)
./run.sh make                        # собрать ./mendrive
./run.sh make test                   # сравнение с исходным ядром + тесты полей
./run.sh make smoke
./run.sh make plot-linear            # эталон TE + сравнение с точным решением
./run.sh make plot CFG=examples/fields_export_preisach_magnetic.cfg OUT=out/p_mag
./run.sh                             # интерактивная оболочка
```
Каталог `MenDriveCpp/` монтируется в контейнер целиком, поэтому `make reflib` видит исходный `../mendrive_core.cpp`, а результаты (`out/`, PNG, CSV) остаются на хосте. Файлы создаются с вашим UID.

Долгий прогон в фоне:
```sh
docker run -d --name p2d -u $(id -u):$(id -g) -v "$(cd ../.. && pwd)":/work/MenDriveCpp \
  -w /work/MenDriveCpp/standalone mendrive-standalone \
  ./mendrive -c examples/preisach2d_magnetic.cfg n_grid=201 threads=6
docker logs -f p2d
```
Прерванный прогон продолжается с `--resume`.

Альтернатива — `docker compose run --rm mendrive make test`.
