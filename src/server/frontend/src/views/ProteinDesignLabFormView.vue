<script>

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
        },
    }
}
</script>

<template>
    <div class="accordion" id="accordionExample">
        <div class="accordion-item">
            <h2 class="accordion-header justify-content-center">
                <button class="accordion-button center-block" :class="{ collapsed: !parametersVisible }"
                    type="button" @click="parametersVisibleShow">
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
                                    data-toggle="tooltip" data-placement="top" data-delay="0" title="Confirm your chains have the residue numbers you're looking to diffuse over. 100 - Diffuses a monomer 100 residues long. 50-100 - Diffuses a hetero-oligomer of lengths 50 and 100. 5-15/A10-25/30-40 - Builds 5-15 residues N-terminally of A10-25 from the input pdb, followed by 30-40 residues to its C-terminus. B1-100/0 100-100 - Generates 100 residue long binders to residues 1-100 of chain B."
                                        id="contigInput" name="contigInput"
                                        aria-describedby="contigInputHelp" placeholder="Contig">
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="timestepsInput" class="col-sm-4 col-form-label">Timesteps</label>
                                <div class="col-sm-8">
                                    <input type="number" required value="50" min="1" step="1" class="form-control"
                                    data-toggle="tooltip" data-placement="top" data-delay="0" title="Desired iterations to generate structure."
                                        id="timestepsInput" name="timestepsInput"
                                        aria-describedby="timestepsInputHelp" placeholder="Timesteps">
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="hotspotsInput" class="col-sm-4 col-form-label">Hotspots (optional)</label>
                                <div class="col-sm-8">
                                    <input type="text" min="0" step="0.001" class="form-control"
                                        id="hotspotsInput" name="hotspotsInput" aria-describedby="hotspotsInputHelp"
                                        data-toggle="tooltip" data-placement="top" data-delay="0" title="The model optionally readily learns that it should be making an interface which involving these hotspot residues. Input is ChainResidueNumber: A100 for residue 100 on chain A."
                                        placeholder="A5,A6,A7..." value="">
                                </div>
                            </div>
                            <div class="mt-1 form-group row">
                                <label for="numOfDesignsInput" class="col-sm-4 col-form-label">Number of designs</label>
                                <div class="col-sm-8">
                                    <input type="number" min="1" step="1" class="form-control"
                                        id="numOfDesignsInput" name="numOfDesignsInput" aria-describedby="numOfDesignsInputHelp"
                                        placeholder="Number of designs" value="1">
                                </div>
                            </div>
                        </div>
                        <div class="justify-content-center col-md-5"
                            style="max-height: 100%;height:100%; overflow-y:scroll">
                            <label for="proteinFileInput" class="col-form-label fs-4">Binder / Scaffolding Structure</label>
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