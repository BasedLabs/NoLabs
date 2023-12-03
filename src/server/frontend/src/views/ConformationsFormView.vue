<script>
export default {
    props: ['onFormSubmit']
}
</script>

<template>
    <form enctype="multipart/form-data" id="inferenceInputFormConformations" class="row" v-on:submit.prevent="onFormSubmit">
        <div class="justify-content-center col-md-5">
            <label for="proteinFileInput" class="col-form-label fs-4">.pdb file with molecule conformations</label>
            <input type="file" required multiple="multiple" accept="text/x-pdb" class="form-control" id="proteinFileInput"
                name="proteinFileInput">
        </div>
        <div class="justify-content-center col-md-7 mt-2">
            <h4>Simulation parameters</h4>
            <div class="form-group row">
                <label for="totalFramesInput" class="col-sm-4 col-form-label">Total frames</label>
                <div class="col-sm-8">
                    <input type="number" required min="0" step="1" class="form-control" id="totalFramesInput"
                    name="totalFramesInput"
                        aria-describedby="totalFramesInputHelp" placeholder="Total frames" value="20000">
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="takeFrameEveryInput" class="col-sm-4 col-form-label">Take frame every</label>
                <div class="col-sm-8">
                    <input type="number" required min="0" step="1" class="form-control" id="takeFrameEveryInput"
                    name="takeFrameEveryInput"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Total frames" value="1000">
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="integratorSelect" class="col-sm-4 col-form-label">Integrator</label>
                <div class="col-sm-8">
                    <select class="form-select" required aria-label="Default select example" id="integratorSelect"
                    name="integratorSelect">
                        <option selected value="LangevinIntegator">Langevin integator</option>
                        <option value="LangevinMiddleIntegrator">Langevin middle integator</option>
                        <option value="NoseHooverIntegrator">Nose hoover integrator</option>
                        <option value="BrownianIntegrator">Brownian integrator</option>
                        <option value="VariableVerletIntegrator">Variable verlet integrator</option>
                    </select>
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="frictionCoeffInput" class="col-sm-4 col-form-label">Friction coeff / 10<sup>-12</sup></label>
                <div class="col-sm-8">
                    <input type="number" min="0" required step="0.01" class="form-control" id="frictionCoeffInput"
                    name="frictionCoeffInput"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Friction coeff" value="1">
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="stepSizeInput" class="col-sm-4 col-form-label">Step size / 10<sup>-12</sup></label>
                <div class="col-sm-8">
                    <input type="number" min="0" required step="0.001" class="form-control" id="stepSizeInput"
                    name="stepSizeInput"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" value="0.002">
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="tempInput" class="col-sm-4 col-form-label">System temp in K</label>
                <div class="col-sm-8">
                    <input type="number" min="0" required step="0.001" class="form-control" id="tempInput"
                    name="tempInput"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" value="273.15">
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="replaceNonstandardResiduesInput" class="col-sm-4 col-form-label">Replace nonstandard residues</label>
                <div class="col-sm-8">
                    <input type="checkbox" class="form-check-input" id="replaceNonstandardResiduesInput"
                    name="replaceNonstandardResiduesInput"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" checked>
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="addMissingAtomsCheckbox" class="col-sm-4 col-form-label">Add missing atoms</label>
                <div class="col-sm-8">
                    <input type="checkbox" class="form-check-input" id="addMissingAtomsCheckbox"
                    name="addMissingAtomsCheckbox"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" checked>
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="addMissingHydrogensCheckbox" class="col-sm-4 col-form-label">Add missing hydrogens</label>
                <div class="col-sm-8">
                    <input type="checkbox" class="form-check-input" id="addMissingHydrogensCheckbox"
                    name="addMissingHydrogensCheckbox"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" checked>
                </div>
            </div>
            <div class="mt-1 form-group row">
                <label for="ignoreMissingAtomsCheckbox" class="col-sm-4 col-form-label">Ignore missing atoms</label>
                <div class="col-sm-8">
                    <input type="checkbox" class="form-check-input" id="ignoreMissingAtomsCheckbox"
                    name="ignoreMissingAtomsCheckbox"
                        aria-describedby="takeFrameEveryInputHelp" placeholder="Simulation step size" checked>
                </div>
            </div>
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
<hr/>

    </form>
</template>