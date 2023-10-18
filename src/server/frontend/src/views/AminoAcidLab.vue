<script>
import factory from '../services/experimentsViewDataFactory';
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import AminoAcidLabFormView from './AminoAcidLabFormView.vue';
import AminoAcidLabInferenceView from './AminoAcidLabInferenceView.vue';

export default {
    data() {
        return {
            experimentsState: store.state.aminoAcid,
            experimentsApi: api.aminoAcid
        }
    },
    components: {
        ExperimentsView
    }
}
</script>

<template>
    <ExperimentsView experimentsState="experimentsState" experimentsApi="experimentsApi">
        <template v-slot:labTitle>
            <h4>Protein lab</h4>
        </template>
        <template v-slot:labForm>
            <AminoAcidLabFormView/>
        </template>
        <template v-slot:lab>
            <AminoAcidLabInferenceView/>
        </template>
    </ExperimentsView >
</template>