# Astro Plotting Server

For use with @shengjiex98's [skynet plotting tool](https://github.com/shengjiex98/skynet-plotting).

## Usage

1. Create a virtual environment and activate it
```bash
$ python --version
3.10.1
$ python -m venv --upgrade-deps venv
$ source ./venv/bin/activate
```

2. Install dependencies
```bash
(venv) $ pip install -r requirements.txt
```

3. Run
```bash
(venv) $ python main.py
```

4. Try it! Example URL: <http://localhost:5000/data?age=6.80&metallicity=-0.35&filters=[%22uprime%22,%22H%22,%22J%22]>
