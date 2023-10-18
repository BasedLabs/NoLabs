<script setup>
import store from '../../storage.js';

const solubility = {
    value: `This protein is soluble with ${+(Math.round(store.state.aminoAcid.inference.solubility * 100.0 + "e+2") + "e-2")}% probability`
};
</script>

<template>
    <div class="col-md-12 row mt-5" id="solubility">
        <h1 id="solubilityText">{{ solubility.value }}</h1>
    </div>
</template>
