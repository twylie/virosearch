import sys


###############################################################################
#                                StdoutMessage                                #
###############################################################################

# A general class for to display messages to the user from STDOUT. This is a
# controlled way to make messages uniform and to also provide a method for a
# fatal message.


class StdoutMessage:

    def __init__(
            self,
            command='',
            function='',
            message='',
            fatal=False,
            verbose=True
    ):

        self.command = command
        self.function = function
        self.message = message
        self.fatal = fatal
        self.verbose = verbose

        message_fields = [
            f'[virosearch {self.command}]',
            f'({self.function})',
            f'--- {self.message}'
        ]

        stdout = ' '.join(message_fields)

        # We can turn off all printing by setting verbose=False.

        if self.verbose is True:
            if self.fatal is True:
                print(stdout)
                print('(EXITING...)')
                sys.exit()
            else:
                print(stdout)

        return


# __END__
