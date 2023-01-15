import time
# import colorama


def log_information(message):
    t = time.strftime('%H:%M:%S')
    print(f'[\033[32m{t}\033[0m] {message}')


def log_warning(message):
    t = time.strftime('%H:%M:%S')
    print(f'[\033[33m{t}\033[0m] {message}')


def log_error(message):
    t = time.strftime('%H:%M:%S')
    print(f'[\033[31m{t}\033[0m] {message}')
