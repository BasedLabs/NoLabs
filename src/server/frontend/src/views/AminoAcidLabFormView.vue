<script>
export default {
    props: {
        onFormSubmit: undefined // method
    }
}
</script>

<template>
    <form enctype="multipart/form-data" id="inferenceInputFormAminoAcid" v-on:submit.prevent="onFormSubmit">
        <div class="row justify-content-center">
            <div class="col-md-6">
                <label for="inputSequence" class="col-form-label fs-4">Paste amino-acid
                    sequence</label>
                <input type="text" pattern="[SMIREATZLHJXGWKBDPCFNQYV]*" class="form-control"
                    value="AAACGAGGCAA" id="inputSequence" name="inputSequence"
                    title="Sequence is only allowed to have {'S', 'M', 'I', 'R', 'E', 'A', 'T', 'Z', 'L', 'H', 'J', 'X', 'G', 'W', 'K', 'B', 'D', 'P', 'C', 'F', 'N', 'Q', 'Y', 'V'} tokens.">
            </div>
        </div>
        <div class="row justify-content-center">
            <div class="col-md-6">
                <label for="inputSequenceFile" class="col-form-label fs-4">Or paste a .FASTA file
                    with
                    amino-acid
                    sequence</label>
                <input type="file" multiple="multiple" accept="text/x-fasta" class="form-control"
                    id="inputSequenceFile" name="inputSequenceFile">
            </div>
        </div>
        <div class="row justify-content-center m-4">
            <div class="col-md-4">
                <button type="submit" id="submitInference" class="btn btn-lg btn-primary">
                    Submit
                </button>
            </div>
        </div>
    </form>
</template>