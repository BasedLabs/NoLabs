describe('Main page', () => {
    const inputPdbId = '#proteinFileInput';
    const inputAminoAcidSequenceId = '#inputSequence';
    const submitAminoAcidSequenceId = '#submitInference';
    const localisationImageClass = '.localisation-image';
    const foldingSelector = '#viewport div';

    it('Opened main page', () => {
        cy.visit('http://127.0.0.1:5000');
    });

    it('Amino Acid lab opened', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid');
    });

    it('Amino-acid sequence input is visible', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid');
        cy.get(inputAminoAcidSequenceId).should('exist');
    });

    it('Inference result is visible', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid');
        cy.get(inputAminoAcidSequenceId).click();
        cy.get(inputAminoAcidSequenceId).type(`AAACGAGGCAA`);
        cy.get(submitAminoAcidSequenceId).click();
        cy.get(foldingSelector, {timeout: 90000}).should('exist');
    });
})