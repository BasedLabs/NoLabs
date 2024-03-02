import subprocess

from fastapi import FastAPI

app = FastAPI(
    title='Reinvent api'
)

@app.post("/run")
async def run_pdb_fixer_endpoint() -> str:
    program = ["reinvent", "-l", "sampling.log", "reinvent/REINVENT4/configs/toml/sampling.toml"]
    res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    print(res.stderr.decode('utf-8'))
    print(res.stdout.decode('utf-8'))
    return 'hello'
