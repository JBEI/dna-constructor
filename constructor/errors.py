class InterpreterError(Exception):
	def __init__(self, value, protocolId):
		self.value = 'Interpreter Error:\n' + value

		self.protocolId = protocolId
	def __str__(self):
		return repr(self.value)

class ConstructError(Exception):
	def __init__(self, value, protocolId):
		self.value = 'Construction Error:\n' + value + '\nIt\'s probably Nick\'s fault- find him or email him at elsbree@lbl.gov for help.'

		self.protocolId = protocolId
	def __str__(self):
		return repr(self.value)
