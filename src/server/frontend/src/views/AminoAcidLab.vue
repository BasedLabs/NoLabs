<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import AminoAcidLabFormView from './AminoAcidLabFormView.vue';
import AminoAcidLabInferenceView from './AminoAcidLabInferenceView.vue';
import ProteinList from '../components/AminoAcidLab/ProteinList.vue';

export default {
    data() {
        return {
            experimentsState: store.state.aminoAcid,
            experimentsApi: api.aminoAcid
        }
    },
    components: {
        ExperimentsView,
        AminoAcidLabFormView,
        AminoAcidLabInferenceView,
        ProteinList
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h4>Protein lab</h4>
        </template>
        <template v-slot:labForm="labForm">
            <AminoAcidLabFormView :onFormSubmit="labForm.onFormSubmit"/>
        </template>
        <template v-slot:experimentMetaData="experimentMetaData" :api="experimentsApi">
            <ProteinList :experiment="experimentMetaData.experiment" :api="experimentsApi"/>
        </template>
        <template v-slot:lab="lab">
            <AminoAcidLabInferenceView :experiment="lab.experiment"/>
        </template>
    </ExperimentsView >
</template>