Creating binstar token:
```
#!/bin/bash
BINSTAR_TOKEN=`binstar auth --create -o omnia --name pdbfixer-travis -s "api:read api:write"`
travis encrypt -r pandegroup/pdbfixer BINSTAR_TOKEN="$BINSTAR_TOKEN"
```
