class NoFileException(Exception):
    def __init__(self):
        self.message = "File Does Not Exist"
        super().__init__(self.message)
