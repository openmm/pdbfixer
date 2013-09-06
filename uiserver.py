from threading import Thread
from SocketServer import ThreadingMixIn
from BaseHTTPServer import HTTPServer, BaseHTTPRequestHandler
from urlparse import parse_qs
import cgi

class _Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.hasSentResponse = False
        if callback is not None:
            queryStart = self.path.find('?')
            if queryStart > -1:
                parameters = parse_qs(self.path[queryStart+1:])
            else:
                parameters = {}
            result = callback(parameters, self)
        if not self.hasSentResponse:
            self.sendResponse(content)
    
    def do_POST(self):
        self.hasSentResponse = False
        if callback is not None:
            parameters = cgi.FieldStorage(fp=self.rfile, headers=self.headers, environ={'REQUEST_METHOD':'POST', 'CONTENT_TYPE':self.headers['Content-Type']})
            callback(parameters, self)
        if not self.hasSentResponse:
            self.sendResponse(content)
    
    def sendResponse(self, response):
        self.hasSentResponse = True
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.send_header("Content-length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)
    
    def sendDownload(self, download, filename):
        self.hasSentResponse = True
        self.send_response(200)
        self.send_header("Content-type", "text/plain")
        self.send_header("Content-length", str(len(download)))
        self.send_header("Content-Disposition", 'attachment; filename="%s"' % filename)
        self.end_headers()
        self.wfile.write(download)

class _ThreadingHTTPServer(ThreadingMixIn, HTTPServer):
    pass

content = ""
callback = None
server = _ThreadingHTTPServer(("localhost", 8000), _Handler)

def beginServing():
    Thread(target=server.serve_forever).start()

def setContent(newContent):
    global content
    content = newContent

def setCallback(newCallback):
    global callback
    callback = newCallback

