Run docker with Rosettafold USERNAME:PASSWORD as following

```bash
docker buildx build --build-arg ROSETTACOMMONS_CONDA_USERNAME=username --build-arg ROSETTACOMMONS_CONDA_PASSWORD=password --progress plain -f build/Dockerfile -t rosettafold .
```

then 

```bash
docker run --name rosettafold -v /mnt/e/RoseTTAFold/bfd:/RoseTTAFold/bfd -v /mnt/e/RoseTTAFold/UniRef30_2020_06:/RoseTTAFold/UniRef30_2020_06 -v /mnt/e/RoseTTAFold/pdb100_2021Mar03:/RoseTTAFold/pdb100_2021Mar03 --gpus=all --net=host rosettafold --port=8000
```