<script>
export default {
    props: ['onFormSubmit']
}
</script>

<template>
    <form enctype="multipart/form-data" id="inferenceInputFormConformations"
        v-on:submit.prevent="onFormSubmit">
        <div class="row justify-content-center">
            <div class="col-md-6">
                <label for="proteinFileInput" class="col-form-label fs-4">.pdb file with molecule conformations</label>
                <input type="file" multiple="multiple" accept="text/x-pdb" class="form-control" id="proteinFileInput"
                    name="proteinFileInput">
            </div>
        </div>
        <div class="row justify-content-center m-4">
            <div class="col-md-4">
                <div class="btn-group mt-4" role="group"
                    aria-label="Basic checkbox toggle button group">
                    <button type="submit" id="submitInference" class="btn btn-lg btn-primary">
                        Submit
                    </button>
                </div>
            </div>
        </div>
    </form>
</template>