<script lang="ts">
	import axios from "axios";
	export let name: string = "Roman";

	let rand = -1;
	function getRand() {
		fetch("./rand")
			.then((r) => r.text())
			.then((d) => (rand = d));
	}

	let files;

	$: if (files) {
		// Note that `files` is of type `FileList`, not an Array:
		// https://developer.mozilla.org/en-US/docs/Web/API/FileList
		console.log(files);

		for (const file of files) {
			console.log(`${file.name}: ${file.size} bytes`);
		}
	}

	function uploadFile(e) {
		e.preventDefault();
		let file = this.state.fileToBeSent;
		const formData = new FormData();

		formData.append("file", file);

		axios
			.post("/api/upload", formData)
			.then((res) => console.log(res))
			.catch((err) => console.warn(err));
	}

	let foo: string = "foo";
	let bar: string = "bar";
	let result = null;
	async function doPost() {
		const res = await fetch("./post", {
			method: "POST",
			body: JSON.stringify({
				foo,
				bar,
			}),
		});
		console.log(res);
		const json = await res.json();
		console.log(json);
		result = JSON.stringify(json);
	}
</script>

<main>
	<div>
		<h1>Hello {name}!</h1>
		<p>
			Visit the <a href="https://svelte.dev/tutorial">Svelte tutorial</a> to
			learn how to build Svelte apps.
		</p>
	</div>

	<div>
		<p>
			Your number is {rand}!
			<button on:click={getRand}>Get a random number</button>
		</p>
	</div>

	<div id="upload">
		<label for="upload_lower"> Upload a structure: </label>
		<input
			accept=".xyz, .cif, .in"
			bind:files
			id="upload_lower"
			name="lower"
			type="file"
		/>
		<button on:click={uploadFile}> Upload </button>
	</div>

	<input bind:value={foo} />
	<input bind:value={bar} />
	<button type="button" on:click={doPost}> Post it. </button>
	<p>Result:</p>
	<pre>
	{result}
	</pre>
</main>

{#if files}
	<h2>Selected files:</h2>
	{#each Array.from(files) as file}
		<p>{file.name} ({file.size} bytes)</p>
	{/each}
{/if}

<style>
	main {
		text-align: center;
		padding: 1em;
		max-width: 240px;
		margin: 0 auto;
	}

	#upload_lower {
		display: flex;
		align-items: center;
		justify-content: center;
		flex-flow: column;
	}

	h1 {
		color: #ff3e00;
		text-transform: uppercase;
		font-size: 4em;
		font-weight: 100;
	}

	@media (min-width: 640px) {
		main {
			max-width: none;
		}
	}
</style>
