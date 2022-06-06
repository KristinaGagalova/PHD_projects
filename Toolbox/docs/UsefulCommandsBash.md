## More tricks and solutions in bash

### Debugging 
```
#start
set -x
#stop
set +x
```

### Weird characters in text file

If you don't have dos2unix installed use this command to remove unuasual characters (\r)

```
sed -i 's/\r$//' filename
```

### Grep tricks

Make it faster when grepping from 2 files
```
grep -wF -f myfile1.txt myfile2
```

### Awk tricks


