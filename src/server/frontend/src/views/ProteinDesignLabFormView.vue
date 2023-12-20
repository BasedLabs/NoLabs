<script>
import socket from '../sockets';

export default {
    props: ['onFormSubmit'],
    data: () => {
        return {
            parametersVisible: true
        }
    },
    methods: {
        parametersVisibleShow() {
            this.parametersVisible = !this.parametersVisible;
        }
    }
}
</script>

<template>
    <div class="accordion" id="accordionExample">
        <div class="accordion-item">
            <h2 class="accordion-header justify-content-center">
                <button class="accordion-button center-block" :class="{ collapsed: !parametersVisible }"
                    type="button" @click="parametersVisible">
                    Protein design parameters
                </button>
            </h2>
            <div class="accordion-collapse collapse" :class="{ show: parametersVisible }">
                <div class="accordion-body">
                    <form enctype="multipart/form-data" id="inferenceInputForm" class="row"
                        v-on:submit.prevent="onFormSubmit">
                        <div class="justify-content-center col-md-7 mt-2">
                            <div class="form-group row">
                                <label for="contigInput" class="col-sm-4 col-form-label">Contig</label>
                                <div class="col-sm-8">
                                    <input type="text" required min="0" step="1" class="form-control"
                                        id="contigInput" name="contigInput"
                                        aria-describedby="contigInputHelp" placeholder="Contig" value="20000">
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="symmetrySelect" class="col-sm-4 col-form-label">Integrator</label>
                                <div class="col-sm-8">
                                    <select class="form-select" required aria-label="Default select example"
                                        id="symmetrySelect" name="symmetrySelect">
                                        <option selected value="">None</option>
                                        <option selected value="cyclic">Cyclic</option>
                                        <option value="dihedral">Dihedral</option>
                                        <option value="tetrahedral">Tetrahedral</option>
                                        <option value="octahedral">Octahedral</option>
                                        <option value="icosahedral">Icosahedral</option>
                                    </select>
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="timestepsInput" class="col-sm-4 col-form-label">Timesteps</label>
                                <div class="col-sm-8">
                                    <input type="number" min="0" required step="1" class="form-control"
                                        id="timestepsInput" name="timestepsInput"
                                        aria-describedby="timestepsInputHelp" placeholder="Friction coeff" value="1">
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="hotspotsInput" class="col-sm-4 col-form-label">Hotspots</label>
                                <div class="col-sm-8">
                                    <input type="text" min="0" required step="0.001" class="form-control"
                                        id="hotspotsInput" name="hotspotsInput" aria-describedby="hotspotsInputHelp"
                                        placeholder="Hotspots" value="">
                                </div>
                            </div>
                        </div>
                        <div class="justify-content-center col-md-5"
                            style="max-height: 100%;height:100%; overflow-y:scroll">
                            <label for="proteinFileInput" class="col-form-label fs-4">Source protein (.pdb)</label>
                            <input type="file" required multiple="multiple" accept="text/x-pdb" class="form-control"
                                id="proteinFileInput" name="proteinFileInput">
                        </div>
                        <div class="row justify-content-center m-4">
                            <div class="col-md-4">
                                <div class="btn-group mt-4" role="group" aria-label="Basic checkbox toggle button group">
                                    <button type="submit" id="submitInference" class="btn btn-lg btn-primary">
                                        Submit
                                    </button>
                                </div>
                            </div>
                        </div>
                    </form>
                </div>
            </div>
        </div>
    </div>
</template>