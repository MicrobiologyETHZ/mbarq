import subprocess


def check_call(command: object) -> int:
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command:
    :return:
    """
    returncode = -1
    try:
        returncode: int = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
    return returncode
