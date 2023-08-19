import time


def header():
    print(time.strftime("%b %d %Y %H:%M:%S", time.localtime()))
    print('')


def print_header(func):
    def wrapper(*args, **kwargs):
        header()
        return func(*args, **kwargs)

    return wrapper


def print_run_info(info):
    def wrapper(func):
        def _wrapper(*args, **kargs):
            print(info + '! Running ...\n')
            res = func(*args, **kargs)
            print(info + ' OK!')
            print('')
            return res

        return _wrapper

    return wrapper


if __name__ == '__main__':
    header()
