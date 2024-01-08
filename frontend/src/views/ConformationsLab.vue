<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import ConformationsFormView from './ConformationsFormView.vue';
import ConformationsInferenceView from './ConformationsInferenceView.vue';

export default {
    data() {
        return {
            experimentsState: store.state.conformations,
            experimentsApi: api.conformations
        }
    },
    components: {
        ExperimentsView,
        ConformationsFormView,
        ConformationsInferenceView
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h2>Conformations lab (pre-alpha)</h2>
        </template>
        <template v-slot:labForm="labForm">
            <ConformationsFormView :onFormSubmit="labForm.onFormSubmit"/>
        </template>
        <template v-slot:lab="lab">
            <ConformationsInferenceView :experiment="lab.experiment"/>
        </template>
    </ExperimentsView >
</template>