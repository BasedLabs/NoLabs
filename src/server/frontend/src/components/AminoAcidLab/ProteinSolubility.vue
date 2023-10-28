<script>
export default {
    props: ['experiment'],
    data() {
        const solubility = `This protein is soluble with ${+(Math.round(this.experiment.data.solubility * 100.0 + "e+2") + "e-2")}% probability`

        return {
            solubility
        }
    }
}
</script>

<template>
    <div class="col-md-12 row mt-5" id="solubility">
        <h1 id="solubilityText">{{ solubility }}</h1>
    </div>
</template>
