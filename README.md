# Требования

- dotnet 9
- Драйвер OpenCL (Nvidia поставляет его вместе с Cuda)

# Клонирование

Нужно склонировать этот репозиторий вместе с его модулями, например:

`git clone --recurse-submodules https://github.com/RocketRide9/course_opencl.git`

# Запуск

`dotnet run` не работает, потому что Visual Studio (Code) запускает проект немного по-другому. Вместо этого нужно:

## Windows

- В Visual Studio (Code) нажать на кнопку запуска

ИЛИ

- В программе "git bash" запустить скрипт: `./Course/run.sh --release`.

## Linux

- В Visual Studio Code тоже сработает по нажатию кнопки.

ИЛИ

- В терминале: `./Course/run.sh --release`.

# Результат тестирований

В `Course/measurements` появится txt файл.
