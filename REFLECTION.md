# Reflection and steps taken to solve the issue.

1. Explore the assignment

2. Explore the repo

Saw .python-version with 3.12.11, so I created an uv venv with it, to prevent any compatibility issues:

```
C:\Users\piete\Repos\pkuppens\enpicom-quality-per-position-challenge>py -V
Python 3.13.5

C:\Users\piete\Repos\pkuppens\enpicom-quality-per-position-challenge>uv venv --python 3.12.11
...

C:\Users\piete\Repos\pkuppens\enpicom-quality-per-position-challenge>.venv\Scripts\activate

(enpicom-quality-per-position-challenge) C:\Users\piete\Repos\pkuppens\enpicom-quality-per-position-challenge>py -V
Python 3.12.11
```

3. start on the assignment

We need the parser first, but from specification, we need to support but text and gzipped text.
Create a function to handle both cases and return a single return type: TextIO

4. Then read per 4 lines, better LINES_PER_READ as this is specified

5. Extract ID, quality string - convert quality string to bytes,
because we'll use a histogram for efficiency