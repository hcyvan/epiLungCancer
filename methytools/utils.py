def is_gz_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            magic_number = file.read(2)
            return magic_number == b'\x1f\x8b'
    except IOError:
        return False